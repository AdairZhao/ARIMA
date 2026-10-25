import Link from "next/link";

export function SiteFooter() {
  return (
    <footer className="site-footer">
      <div className="site-footer__inner">
        <div>
          <p className="brand brand--footer">
            <span className="brand__mark">東</span>
            <span className="brand__name">EASTBOUND</span>
          </p>
          <p className="muted">
            Practical China travel for first-time visitors — visas, payments,
            trains, and routes that actually work on the ground.
          </p>
        </div>
        <div className="footer-links">
          <Link href="/plan">Plan</Link>
          <Link href="/destinations">Destinations</Link>
          <Link href="/itineraries">Itineraries</Link>
          <Link href="/transport">Transport</Link>
          <Link href="/essentials">Essentials</Link>
        </div>
      </div>
      <p className="site-footer__note">
        English-first guide · Built for travelers planning from abroad
      </p>
    </footer>
  );
}
